% ─── Physical constants ──────────────────────────────────────────────────────
G         = 9.80665;          % m/s²
P_ATM     = 101325.0;         % Pa
RHO_WATER = 1000.0;           % kg/m³
RHO_AIR   = 1.200;            % kg/m³  (Rzeszow site, FDR §4.8.3)
GAMMA     = 1.4;              % air adiabatic exponent
CD_NOZZLE = 0.82;             % nozzle discharge coefficient (FDR §4.2.2)

% ─── Vehicle — Ad Astra Gades MkIV  (masses from FDR §4.1.3) ────────────────
AIRFRAME     = 0.891 + 0.444;
MASS_PAYLOAD = 0.800;
DRY_MASS     = AIRFRAME + MASS_PAYLOAD;

N_STAGES    = 2;
V_BOTTLE    = 0.002;
WATER_STAGE = MASS_PAYLOAD / (N_STAGES * RHO_WATER);
P_FILL      = 6.0 * P_ATM;
SHELL_MASS  = 0.070;

NOZZLE_D  = 0.009;
A_NOZZLE  = pi * (NOZZLE_D / 2)^2;

% Launch tube slides up inside the bottle through the nozzle bore, so its
% outer diameter must be smaller than the nozzle exit — the piston area is the
% TUBE cross-section, not the wide PCO neck.
PCO_NECK_D = 0.0217;          % m — PCO neck bore (reference only)
TUBE_OD    = 0.008;           % m — tube outer diameter (< NOZZLE_D = 9 mm)
TUBE_LEN   = 0.25;
A_TUBE     = pi * (TUBE_OD / 2)^2;

% ─── Thrust calibration ──────────────────────────────────────────────────────
% A lumped efficiency factor on the water thrust that folds in the real-world
% impulse the simple incompressible-jet model misses (air-pulse tail, etc.).
% Burnout speed ≈ 8.8*THRUST_CAL − 5.2 m/s, so:
%   1.5 -> ~8 m/s -> ~5 m apogee | 2.0 -> ~12.5 m/s -> ~10 m | 2.2 -> ~12 m
THRUST_CAL = 2.00;

BODY_D  = 0.105;
A_BODY  = pi * (BODY_D / 2)^2;
CD_BODY = 0.50;

INCLINATION = deg2rad(85.0);
cos_el = cos(INCLINATION);
sin_el = sin(INCLINATION);
V_RAIL = 2.0;                 % m/s — rail holds launch attitude below this speed

% ─── Derived air volumes ─────────────────────────────────────────────────────
AIR_FULL = V_BOTTLE - WATER_STAGE;
AIR_TUBE = AIR_FULL - A_TUBE * TUBE_LEN;

M0 = DRY_MASS + N_STAGES * (WATER_STAGE * RHO_WATER + SHELL_MASS);

% Bundle constants into structs so local functions stay clean
C.G = G; C.P_ATM = P_ATM; C.RHO_WATER = RHO_WATER; C.RHO_AIR = RHO_AIR;
C.GAMMA = GAMMA; C.CD_NOZZLE = CD_NOZZLE; C.V_BOTTLE = V_BOTTLE;
C.THRUST_CAL = THRUST_CAL;
C.A_NOZZLE = A_NOZZLE; C.A_TUBE = A_TUBE; C.P_FILL = P_FILL;
C.AIR_TUBE = AIR_TUBE; C.AIR_FULL = AIR_FULL;
C.CD_BODY = CD_BODY; C.A_BODY = A_BODY;
C.cos_el = cos_el; C.sin_el = sin_el;
C.V_RAIL = V_RAIL;
C.M0 = M0;

% ─── Phase 0 — launch tube ───────────────────────────────────────────────────
% State: [d; v; air].  The thin tube (OD < nozzle) only boosts the rocket if its
% piston force beats gravity on the rail; for this heavy glider it does not, so
% the phase is skipped and the bottle launches from rest under water thrust.
a_tube0 = (P_FILL - P_ATM) * A_TUBE / M0 - G * sin_el;   % peak tube accel
if a_tube0 > 0.0
    opts_tube = odeset('Events',  @(t,y) ev_tube_exit(t, y, TUBE_LEN), ...
                       'MaxStep', 5e-5, 'RelTol', 1e-7, 'AbsTol', 1e-9);
    [t_ube, y_ube] = ode45(@(t,y) tube_ode(t, y, C), ...
                            [0.0, 2.0], [0.0; 0.0; AIR_TUBE], opts_tube);
    t_ube       = t_ube';
    d_ube_arr   = y_ube(:,1)';
    v_ube_arr   = y_ube(:,2)';
    air_ube_arr = y_ube(:,3)';
else
    t_ube = 0.0;  d_ube_arr = 0.0;  v_ube_arr = 0.0;  air_ube_arr = AIR_TUBE;
end
v_ube_exit  = v_ube_arr(end);
t_ube_end   = t_ube(end);

% ─── Phase 1 — stage-1 water burn ────────────────────────────────────────────
% State: [x; z; vx; vz; water; mass]
opts_burn = odeset('Events',  @ev_water_out, ...
                   'MaxStep', 5e-5, 'RelTol', 1e-7, 'AbsTol', 1e-9);

x0   = d_ube_arr(end) * cos_el;
z0   = d_ube_arr(end) * sin_el;
y1_0 = [x0; z0; v_ube_exit*cos_el; v_ube_exit*sin_el; WATER_STAGE; M0];

[t1, y1] = ode45(@(t,y) burn_ode(t, y, C.AIR_TUBE, C), ...
                  [t_ube_end, t_ube_end+5.0], y1_0, opts_burn);
t1 = t1';  y1 = y1';
t_burn1_end = t1(end);
m_post1     = y1(6,end) - SHELL_MASS;     % jettison empty stage-1 shell

% ─── Phase 2 — stage-2 water burn ────────────────────────────────────────────
y2_0 = [y1(1:4,end); WATER_STAGE; m_post1];

[t2, y2] = ode45(@(t,y) burn_ode(t, y, C.AIR_FULL, C), ...
                  [t_burn1_end, t_burn1_end+5.0], y2_0, opts_burn);
t2 = t2';  y2 = y2';
t_burn2_end = t2(end);
m_coast     = y2(6,end) - SHELL_MASS;     % jettison empty stage-2 shell

% ─── Phase 3 — ballistic coast to apogee ─────────────────────────────────────
opts_coast = odeset('Events',  @ev_apogee, ...
                    'MaxStep', 0.01, 'RelTol', 1e-7, 'AbsTol', 1e-9);

y3_0 = y2(1:4,end);
[t3, y3] = ode45(@(t,y) coast_ode(t, y, m_coast, C), ...
                  [t_burn2_end, t_burn2_end+60.0], y3_0, opts_coast);
t3 = t3';  y3 = y3';

apogee_z = y3(2,end);
apogee_t = t3(end);

% ─── Reconstruct thrust and mass arrays for plotting ─────────────────────────

% Tube phase
P_ph0    = P_FILL * (AIR_TUBE ./ max(air_ube_arr, 1e-12)).^GAMMA;
T_ph0    = max(P_ph0 - P_ATM, 0.0) .* A_TUBE;
spd_ph0  = v_ube_arr;
z_ph0    = d_ube_arr .* sin_el;
mass_ph0 = M0 * ones(size(t_ube));

% Stage-1 water burn
z1   = y1(2,:);  vx1 = y1(3,:);  vz1 = y1(4,:);
w1   = y1(5,:);  m1  = y1(6,:);
dP1  = max(P_FILL * (AIR_TUBE ./ max(V_BOTTLE - max(w1,0), 1e-12)).^GAMMA - P_ATM, 0.0);
T1   = (w1 > 1e-10) .* (THRUST_CAL * 2 * CD_NOZZLE * A_NOZZLE .* dP1);
spd1 = hypot(vx1, vz1);

% Stage-2 water burn
z2   = y2(2,:);  vx2 = y2(3,:);  vz2 = y2(4,:);
w2   = y2(5,:);  m2  = y2(6,:);
dP2  = max(P_FILL * (AIR_FULL ./ max(V_BOTTLE - max(w2,0), 1e-12)).^GAMMA - P_ATM, 0.0);
T2   = (w2 > 1e-10) .* (THRUST_CAL * 2 * CD_NOZZLE * A_NOZZLE .* dP2);
spd2 = hypot(vx2, vz2);

% Coast
z3   = y3(2,:);  vx3 = y3(3,:);  vz3 = y3(4,:);
spd3 = hypot(vx3, vz3);
T3   = zeros(size(t3));
m3   = m_coast * ones(size(t3));

% Concatenate (tube → stage-1 burn → stage-2 burn → coast)
t_all   = [t_ube,    t1,   t2,   t3  ];
z_all   = [z_ph0,    z1,   z2,   z3  ];
spd_all = [spd_ph0,  spd1, spd2, spd3];
T_all   = [T_ph0,    T1,   T2,   T3  ];
m_all   = [mass_ph0, m1,   m2,   m3  ];

% ─── Console summary ─────────────────────────────────────────────────────────
prop_mass = N_STAGES * WATER_STAGE * RHO_WATER;
fprintf('Liftoff mass       : %.3f kg\n', M0);
fprintf('Water per stage    : %.1f L  (%d%% fill)\n', ...
        WATER_STAGE*1e3, round(WATER_STAGE/V_BOTTLE*100));
fprintf('Tube-exit speed    : %.2f m/s   at t = %.4f s\n', v_ube_exit, t_ube_end);
fprintf('Stage-1 burnout    : %.2f m/s   at t = %.4f s   z = %.2f m\n', ...
        spd1(end), t_burn1_end, z1(end));
fprintf('Stage-2 burnout    : %.2f m/s   at t = %.4f s   z = %.2f m\n', ...
        spd2(end), t_burn2_end, z2(end));
fprintf('Thrust calibration : x%.2f\n', THRUST_CAL);
fprintf('Apogee             : %.2f m        at t = %.2f s\n', apogee_z, apogee_t);
fprintf('Mass at apogee     : %.3f kg  (expected DRY_MASS = %.3f kg)\n', ...
        m_coast, DRY_MASS);

% ─── Regulation compliance ───────────────────────────────────────────────────
verd = {'FAIL', 'PASS'};   % indexed by logical+1
fprintf('\nCompliance (ASRW 2026):\n');
fprintf('  2.1.7 propellant <= payload : %.2f <= %.2f kg   [%s]\n', ...
        prop_mass, MASS_PAYLOAD, verd{(prop_mass <= MASS_PAYLOAD + 1e-9) + 1});
fprintf('  2.1.8 take-off weight <= 5  : %.2f kg            [%s]\n', ...
        M0, verd{(M0 <= 5.0) + 1});
fprintf('  2.2.1 fill pressure <= 10   : %.1f atm           [%s]\n', ...
        P_FILL/P_ATM, verd{(P_FILL <= 10.0*P_ATM) + 1});
fprintf('  Tube OD < nozzle bore       : %.1f < %.1f mm      [%s]\n', ...
        TUBE_OD*1e3, NOZZLE_D*1e3, verd{(TUBE_OD < NOZZLE_D) + 1});

% ─── Plots ───────────────────────────────────────────────────────────────────
phase_t   = [t_ube_end,   t_burn1_end,                 t_burn2_end              ];
phase_lbl = {'Tube exit', 'Stage 1 burnout/jettison',  'Stage 2 burnout/jettison'};
phase_rgb = [0.25 0.41 0.88; 1.00 0.55 0.00; 0.70 0.13 0.13];

figure('Position', [100 100 1100 700]);
sgtitle('Water Rocket Ascent — Ad Astra Gades MkIV  (RK45)', ...
        'FontSize', 13, 'FontWeight', 'bold');

titles    = {'Speed vs Time', 'Altitude vs Time', 'Thrust vs Time', 'Vehicle Mass vs Time'};
ylabels   = {'Speed (m/s)', 'Altitude AGL (m)', 'Thrust (N)', 'Mass (kg)'};
ydata     = {spd_all, z_all, T_all, m_all};
colors    = {[0.27 0.51 0.71], [0.18 0.55 0.34], [0.86 0.08 0.24], [0.60 0.20 0.80]};
dlabels   = {'Speed |v|', ...
             sprintf('Apogee: %.1f m  at t=%.1f s', apogee_z, apogee_t), ...
             'Thrust (N)', 'Total mass'};

for k = 1:4
    ax = subplot(2,2,k);
    plot(ax, t_all, ydata{k}, 'Color', colors{k}, 'LineWidth', 1.8, ...
         'DisplayName', dlabels{k});
    hold(ax, 'on');
    for p = 1:3
        xline(ax, phase_t(p), '--', 'Color', phase_rgb(p,:), 'LineWidth', 1.1, ...
              'Label', phase_lbl{p}, 'LabelHorizontalAlignment', 'right', ...
              'Alpha', 0.8);
    end
    if k == 4
        yline(ax, m_coast, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.1, ...
              'Label', sprintf('Dry mass  %.3f kg', m_coast));
    end
    xlabel(ax, 'Time (s)');  ylabel(ax, ylabels{k});
    title(ax, titles{k});
    legend(ax, 'FontSize', 7);  grid(ax, 'on');
end

% ═══════════════════════════════════════════════════════════════════════════════
% Local functions  (must appear after all script code)
% ═══════════════════════════════════════════════════════════════════════════════

function dydt = tube_ode(~, y, C)
    v   = y(2);
    air = y(3);
    P   = C.P_FILL * (C.AIR_TUBE / max(air, 1e-12))^C.GAMMA;
    dP  = max(P - C.P_ATM, 0.0);
    T   = dP * C.A_TUBE;
    D   = 0.5 * C.RHO_AIR * v^2 * C.CD_BODY * C.A_BODY;
    a   = (T - D) / C.M0 - C.G * C.sin_el;
    dydt = [v; a; C.A_TUBE * max(v, 0.0)];
end

function [val, isterminal, direction] = ev_tube_exit(~, y, tube_len)
    val        = y(1) - tube_len;
    isterminal = 1;
    direction  = 1;
end

function dydt = burn_ode(~, y, air0_ref, C)
    vx = y(3);  vz = y(4);  water = y(5);  mass = y(6);
    v   = hypot(vx, vz);
    air = C.V_BOTTLE - max(water, 0.0);
    P   = C.P_FILL * (air0_ref / max(air, 1e-12))^C.GAMMA;
    dP  = max(P - C.P_ATM, 0.0);
    if dP > 0.0 && water > 1e-10
        v_jet  = sqrt(2.0 * dP / C.RHO_WATER);
        thrust = C.THRUST_CAL * 2.0 * C.CD_NOZZLE * C.A_NOZZLE * dP;
        dwater = -C.CD_NOZZLE * C.A_NOZZLE * v_jet;
        dmass  = C.RHO_WATER * dwater;
    else
        thrust = 0.0;  dwater = 0.0;  dmass = 0.0;
    end
    D = 0.5 * C.RHO_AIR * v^2 * C.CD_BODY * C.A_BODY;
    % Attitude: the guide rail holds the body on the launch angle until the
    % rocket has real airspeed.  Letting thrust follow velocity below ~2 m/s
    % hits the gravity-turn singularity (direction tumbles at near-zero v).
    if v > C.V_RAIL
        ux = vx/v;  uz = vz/v;
    else
        ux = C.cos_el;  uz = C.sin_el;
    end
    ax = (thrust - D) * ux / mass;
    az = (thrust - D) * uz / mass - C.G;
    dydt = [vx; vz; ax; az; dwater; dmass];
end

function [val, isterminal, direction] = ev_water_out(~, y)
    % Fire just above empty: burn_ode freezes dwater below 1e-10, so water can
    % stall at a tiny positive value and never cross exactly zero.
    val        = y(5) - 1e-8;
    isterminal = 1;
    direction  = -1;
end

function dydt = coast_ode(~, y, m_coast, C)
    vx = y(3);  vz = y(4);
    v  = hypot(vx, vz);
    D  = 0.5 * C.RHO_AIR * v^2 * C.CD_BODY * C.A_BODY;
    if v > 1e-6
        ux = vx/v;  uz = vz/v;
    else
        ux = 0.0;  uz = 1.0;
    end
    dydt = [vx; vz; -D*ux/m_coast; -D*uz/m_coast - C.G];
end

function [val, isterminal, direction] = ev_apogee(~, y)
    % Terminate at apogee (vz crosses 0 downward) or, as a safety net, on
    % ground contact — so a bad burnout state can never integrate to -5000 m.
    val        = [y(4); y(2)];
    isterminal = [1; 1];
    direction  = [-1; -1];
end
