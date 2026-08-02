--!strict
-- NavGen.Boundary — outlines from the per-part LOCAL GRIDS.
--
-- WHAT CHANGED, AND WHY IT IS THE WHOLE DESIGN. The previous version of this
-- module read the world-aligned surfel raster from `Floor`. On that raster a
-- rotated part's rim is a jagged run of axis-aligned cells, so the module spent
-- almost all of its length undoing that: a Euclidean distance transform, a
-- greedy total-least-squares line fit to recover the angle the raster had
-- destroyed, corner recovery from where the fit failed, then miter and bevel
-- cleanup. Every one of those stages was a fit, so every one had a tolerance,
-- and the tolerances fought each other.
--
-- `LocalGrid` already solved it. Each block part is sampled on ITS OWN axes, so
-- a rotated part's rim falls on whole lattice lines of its own grid. There is no
-- staircase to fit away, and therefore nothing to fit: the outline of the cell
-- mask IS the outline, exactly, at whatever angle the part is at.
--
-- Two consequences worth stating plainly, because they are what the fitting
-- pipeline could never have:
--
--   * HEIGHTS ARE EXACT BY CONSTRUCTION. A vertex is a point on the part's own
--     face plane, computed as origin + u*uc + v*vc. Nothing is sampled from a
--     neighbouring cell and nothing is interpolated, so a vertex cannot pick up
--     a height from the wrong side of a cliff.
--   * THE RIM IS THE PART'S RIM. The cell lattice covers floor(2*ext/step)
--     whole cells, so the mask stops up to one step short of the real face
--     edge. `LocalGrid` keeps `center`, `uExt` and `vExt` for exactly this
--     reason: an outline vertex sitting on the mask's outermost lattice line is
--     mapped to the part's true extent instead. The polygon lands on the part's
--     face, not on the lattice that approximates it.
--
-- Non-block parts (Unions, MeshParts, wedges) have no meaningful surface axes,
-- so `LocalGrid` gives them a world-aligned fallback grid. Those are traced by
-- the same code in world XZ, and they are the only place a staircased rim can
-- still appear. That is a known limit of the fallback, not of this module.

local Boundary = {}

export type Config = {
	-- Fallback grids only. A world-aligned grid can hold two surfaces at once
	-- (a mesh with a ledge), so its trace still needs to know what counts as a
	-- cliff. Block grids are a single plane and never consult this.
	stepTol: number?,

	-- How far, in cells, a boundary edge may sit from a killer's face plane and
	-- still be considered part of it. A staircased run is quantized by at most
	-- one cell, so slightly over one is the honest window: wide enough to catch
	-- every step of the staircase, too narrow to reach a plane that is not the
	-- one that made this edge.
	snapCells: number?,
}

local DEFAULT = {
	stepTol = 2.2,
	snapCells = 1.1,
}

local function merged(cfg): any
	local c = {}
	for k, v in pairs(DEFAULT) do c[k] = v end
	if cfg then for k, v in pairs(cfg) do if v ~= nil then c[k] = v end end end
	return c
end

--------------------------------------------------------------------------
-- Contour of a cell mask, in whatever integer lattice the mask is indexed by.
--
-- Boundary EDGES are chained rather than boundary cells walked: the boundary of
-- a set of cells is a closed loop by construction, so there is no open end to
-- chase and no tolerance involved.
--
-- `cliff(a, b)` reports that two lattice-adjacent cells are not continuous
-- ground even though both are in the mask. On a block grid it is never true.
--------------------------------------------------------------------------

local DIRS = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } }

-- lattice corners of cell (u,v)'s edge in direction di, oriented so the loop
-- runs with the interior on its left
local function corners(u: number, v: number, di: number)
	if di == 1 then return { u + 1, v }, { u + 1, v + 1 }
	elseif di == 2 then return { u + 1, v + 1 }, { u, v + 1 }
	elseif di == 3 then return { u, v + 1 }, { u, v }
	else return { u, v }, { u + 1, v } end
end

local function traceMask(cells: {any}, index: {[string]: any}, deadIndex: {[string]: any}?, cliff: ((any, any) -> boolean)?)
	local segs: {any} = {}
	local byStart: {[string]: {any}} = {}
	for _, cell in ipairs(cells) do
		for di, d in ipairs(DIRS) do
			local nkey = (cell.ui + d[1]) .. ":" .. (cell.vi + d[2])
			local nb = index[nkey]
			if not nb or (cliff and cliff(cell, nb)) then
				local a, b = corners(cell.ui, cell.vi, di)
				-- WHAT IS ON THE OTHER SIDE OF THIS EDGE, by name. A dead cell was
				-- recorded with the instance that killed it, so the boundary never
				-- has to re-probe the world to ask why the floor stops here -- which
				-- is the whole reason LocalGrid keeps the attribution.
				local dead = deadIndex and deadIndex[nkey]
				local s = { a = a, b = b, cell = cell, killer = dead and dead.killer or nil }
				segs[#segs + 1] = s
				local k = a[1] .. "," .. a[2]
				local bucket = byStart[k]
				if not bucket then bucket = {}; byStart[k] = bucket end
				bucket[#bucket + 1] = s
			end
		end
	end

	-- CHAINING AT A JUNCTION is not a free choice. A cliff between two cells
	-- that are both in the mask emits an edge from EACH side, on the same
	-- lattice edge in opposite directions, because the upper rim and the lower
	-- rim really are two different boundaries that coincide in plan. The mask's
	-- outline therefore has zero-width slits whose ends are lattice vertices
	-- with four boundary edges leaving them. Taking whichever candidate came
	-- first there splices one rim onto the other and drops degenerate slivers.
	--
	-- The rule that makes it deterministic is the standard face traversal of an
	-- embedded planar graph: arriving along an edge, leave on the next edge
	-- CLOCKWISE from it. With the interior kept on the left by `corners`, that
	-- walks each rim whole and turns a slit around at its tip.
	local function angle(s: any): number
		return math.atan2(s.b[2] - s.a[2], s.b[1] - s.a[1])
	end
	local TAU = math.pi * 2
	local used: {[any]: boolean} = {}
	local loops: {{any}} = {}
	for _, s0 in ipairs(segs) do
		if not used[s0] then
			local loop, cur = {}, s0
			while cur and not used[cur] do
				used[cur] = true
				loop[#loop + 1] = cur
				local cands = byStart[cur.b[1] .. "," .. cur.b[2]]
				local nxt, bestTurn = nil, math.huge
				if cands then
					local back = angle(cur) + math.pi
					for _, s in ipairs(cands) do
						if not used[s] then
							local turn = (back - angle(s)) % TAU
							-- straight back is the full turn, so it is taken only
							-- when nothing else is left: the tip of a slit
							if turn <= 1e-9 then turn = TAU end
							if turn < bestTurn then nxt, bestTurn = s, turn end
						end
					end
				end
				cur = nxt
			end
			if #loop >= 4 then loops[#loops + 1] = loop end
		end
	end
	return loops
end

--------------------------------------------------------------------------
-- SNAPPING A RUN TO THE THING THAT CAUSED IT
--
-- The outer rim of a grid is exact because it is the host part's own face. A
-- HOLE is not: it is the footprint of some other part standing on the host, and
-- if that part disagrees with the host's orientation its footprint staircases
-- across the host's lattice exactly the way everything used to staircase across
-- the world's. Measured on SmallMap's ground slab: the two buildings sharing
-- the ground's 8 degree yaw gave holes with 120-stud straight edges, while the
-- two that did not gave holes of 130 and 152 vertices whose longest edge was 2
-- studs and 1 stud. Same map, same code -- the difference is only whether two
-- parts happened to agree.
--
-- Fitting a line to that staircase is what the old world-raster module did, and
-- it is still the wrong answer: the true line is not unknown. `LocalGrid` killed
-- each of those cells and wrote down WHICH INSTANCE killed it. So take that
-- part's side faces, intersect them with the host's face plane, and put the
-- boundary on the resulting line exactly. No fit, no tolerance on the geometry,
-- and the answer is right at any angle because it never involved the lattice.
--
-- A killer's side face, as a line in the host's own 2D face coordinates. A point
-- (a, b) means the world point origin + u*a + v*b, so a face plane through Q
-- with normal m is  a*(u.m) + b*(v.m) = (Q-P0).m  -- an ordinary 2D line.
--------------------------------------------------------------------------

local function sideLines(part: BasePart, origin: Vector3, u: Vector3, v: Vector3): {any}
	local cf, sz = part.CFrame, part.Size
	local axes = {
		{ dir = cf.RightVector, ext = sz.X * 0.5 },
		{ dir = cf.UpVector, ext = sz.Y * 0.5 },
		{ dir = cf.RightVector:Cross(cf.UpVector), ext = sz.Z * 0.5 },
	}
	local out = {}
	for _, ax in ipairs(axes) do
		for _, sgn in ipairs({ 1, -1 }) do
			local m = ax.dir * sgn
			local A, B = u:Dot(m), v:Dot(m)
			-- A face parallel to the host's own plane cuts it in nothing. That is
			-- the killer's top and bottom, and skipping them is what leaves the
			-- four side faces that actually form the wall.
			local mag = math.sqrt(A * A + B * B)
			if mag > 1e-4 then
				local Q = part.Position + m * ax.ext
				out[#out + 1] = { A = A / mag, B = B / mag, C = (Q - origin):Dot(m) / mag }
			end
		end
	end
	return out
end

local function lineDist(L: any, a: number, b: number): number
	return math.abs(L.A * a + L.B * b - L.C)
end

-- Where two of those lines cross. Near-parallel has no usable crossing, and
-- forcing one produces the miter spike that plagued the old module, so it is
-- reported rather than invented.
local function intersect(L1: any, L2: any): (number?, number?)
	local det = L1.A * L2.B - L2.A * L1.B
	if math.abs(det) < 1e-6 then return nil, nil end
	return (L1.C * L2.B - L2.C * L1.B) / det, (L1.A * L2.C - L2.A * L1.C) / det
end

local function signedArea(pts: {Vector3}): number
	local a2 = 0
	for i = 1, #pts do
		local p, q = pts[i], pts[i % #pts + 1]
		a2 += p.X * q.Z - q.X * p.Z
	end
	return a2 * 0.5
end

--------------------------------------------------------------------------
-- One grid -> its rings, in world space
--------------------------------------------------------------------------

local function ringsOfGrid(g: any, c: any, stats: any)
	local isBlock = not g.fallback and g.n ~= nil
	local step = g.step

	-- Lattice index -> host face coordinate, RIM-EXACT at the extremes. The
	-- lattice covers nu whole cells of `step`, which is up to one step short of
	-- the face's real half-extent, so the outermost lattice line is mapped to the
	-- part's own extent instead. This is what LocalGrid keeps center/uExt/vExt
	-- for, and it is the same idea the killer snapping below generalises.
	local toFace
	if isBlock then
		local nu = math.max(1, math.floor(2 * g.uExt / step + 1e-6))
		local nv = math.max(1, math.floor(2 * g.vExt / step + 1e-6))
		toFace = function(iu: number, iv: number): (number, number)
			local uc = (iu <= 0) and 0 or (iu >= nu and 2 * g.uExt or iu * step)
			local vc = (iv <= 0) and 0 or (iv >= nv and 2 * g.vExt or iv * step)
			return uc, vc
		end
	else
		-- A fallback grid has no face rectangle to be exact about: its indices are
		-- world lattice lines, they run negative, and there is no extent to clamp
		-- the outermost one to. Straight through.
		toFace = function(iu: number, iv: number): (number, number)
			return iu * step, iv * step
		end
	end

	local toWorld
	if isBlock then
		toWorld = function(a: number, b: number): Vector3
			return g.origin + g.u * a + g.v * b
		end
	else
		-- Fallback grid: ui/vi are world XZ lattice indices and the surface is not
		-- one plane, so a corner takes the height of the highest cell touching it.
		-- No search and no fit: the four cells at a lattice corner are the only
		-- candidates there are.
		toWorld = function(a: number, b: number): Vector3
			local iu, iv = math.floor(a / step + 0.5), math.floor(b / step + 0.5)
			local bestY = nil
			for du = -1, 0 do
				for dv = -1, 0 do
					local cell = g.index[(iu + du) .. ":" .. (iv + dv)]
					if cell and (not bestY or cell.pos.Y > bestY) then bestY = cell.pos.Y end
				end
			end
			return Vector3.new(a, bestY or 0, b)
		end
	end

	local cliff = nil
	if not isBlock then
		cliff = function(a: any, b: any): boolean
			return math.abs(a.pos.Y - b.pos.Y) > c.stepTol
		end
	end

	-- A killer's side lines, computed once per killer per grid.
	local linesOf: {[Instance]: any} = {}
	local function killerLines(k: Instance?)
		if not k or not isBlock then return nil end
		local cached = linesOf[k]
		if cached ~= nil then return cached ~= false and cached or nil end
		-- Only a Block has faces worth stealing. A Union or a MeshPart has a
		-- bounding box that is not its shape, and snapping a boundary onto a box
		-- that the geometry does not fill would hand back ground that is not
		-- there. Those keep the lattice outline, which is conservative.
		local ok = k:IsA("Part") and (k :: any).Shape == Enum.PartType.Block
		local lines = ok and sideLines(k :: BasePart, g.origin, g.u, g.v) or false
		linesOf[k] = lines
		return lines ~= false and lines or nil
	end

	local snapMax = c.snapCells * step

	local rings: {any} = {}
	for _, loop in ipairs(traceMask(g.cells, g.index, g.deadIndex, cliff)) do
		-- every boundary edge in face coordinates, with the line it belongs to
		local edges = {}
		for _, s in ipairs(loop) do
			local a1, a2 = toFace(s.a[1], s.a[2])
			local b1, b2 = toFace(s.b[1], s.b[2])
			local line = nil
			local cand = killerLines(s.killer)
			if cand then
				-- Only ever the lines of THIS edge's own killer, so a boundary can
				-- never be pulled onto a plane belonging to some unrelated part that
				-- happens to pass nearby. The edge midpoint is the test: on a
				-- staircased run half the edges lie along the true line and half run
				-- across it, and requiring parallelism would throw away the second
				-- half. Both kinds sit within half a cell of the line.
				local mu, mv = (a1 + b1) * 0.5, (a2 + b2) * 0.5
				local best = math.huge
				for _, L in ipairs(cand) do
					local d = lineDist(L, mu, mv)
					if d < best and d <= snapMax then best, line = d, L end
				end
			end
			edges[#edges + 1] = { a = { a1, a2 }, b = { b1, b2 }, line = line }
		end
		if #edges < 4 then continue end

		-- Group consecutive edges sharing a line into runs. Rotate first so a run
		-- never straddles the arbitrary start of the loop.
		local n = #edges
		local start = 1
		for i = 1, n do
			if edges[i].line ~= edges[(i - 2) % n + 1].line then start = i; break end
		end
		local runs = {}
		for k = 0, n - 1 do
			local e = edges[(start - 1 + k) % n + 1]
			local last = runs[#runs]
			if last and last.line and last.line == e.line then
				last.edges[#last.edges + 1] = e
			else
				runs[#runs + 1] = { line = e.line, edges = { e } }
			end
		end

		-- Corners are where two runs meet: two snapped lines cross exactly, and
		-- anything else keeps the lattice corner it already had. A run on a line
		-- contributes only that corner; a run that was not snapped also
		-- contributes its own interior points.
		local face = {}
		for i, run in ipairs(runs) do
			local prev = runs[(i - 2) % #runs + 1]
			local shared = run.edges[1].a
			local pu, pv = shared[1], shared[2]
			if run.line and prev.line and run.line ~= prev.line then
				local ix, iy = intersect(prev.line, run.line)
				if ix and math.abs(ix - pu) <= snapMax * 4 and math.abs(iy - pv) <= snapMax * 4 then
					pu, pv = ix, iy
					stats.corners += 1
				else
					stats.unstableCorners += 1
				end
			elseif run.line then
				-- entering a snapped run from an unsnapped one: slide the lattice
				-- corner onto the line rather than leaving a step at the join
				local d = run.line.A * pu + run.line.B * pv - run.line.C
				pu, pv = pu - run.line.A * d, pv - run.line.B * d
			elseif prev.line then
				local d = prev.line.A * pu + prev.line.B * pv - prev.line.C
				pu, pv = pu - prev.line.A * d, pv - prev.line.B * d
			end
			face[#face + 1] = { pu, pv }
			if not run.line then
				for k = 2, #run.edges do
					local p = run.edges[k].a
					face[#face + 1] = { p[1], p[2] }
				end
			end
		end

		-- drop points that carry no turn (a straight lattice run, or two runs that
		-- resolved onto the same line)
		local kept = {}
		for i = 1, #face do
			local p = face[i]
			local a = face[(i - 2) % #face + 1]
			local b = face[i % #face + 1]
			local d1u, d1v = p[1] - a[1], p[2] - a[2]
			local d2u, d2v = b[1] - p[1], b[2] - p[2]
			if math.abs(d1u * d2v - d1v * d2u) > 1e-6 then kept[#kept + 1] = p end
		end
		if #kept < 3 then kept = face end
		if #kept < 3 then continue end

		local verts = table.create(#kept)
		for i, p in ipairs(kept) do verts[i] = toWorld(p[1], p[2]) end
		rings[#rings + 1] = { verts = verts, area = signedArea(verts) }
	end
	return rings
end

--------------------------------------------------------------------------
-- Boundary.fromLocal — the entry point
--------------------------------------------------------------------------

function Boundary.fromLocal(localData: any, cfg: Config?)
	local c = merged(cfg)
	c.stepTol = c.stepTol or 2.2
	local t0 = os.clock()

	local regions: {any} = {}
	local stats = {
		parts = 0, block = 0, fallback = 0,
		regions = 0, holes = 0, verts = 0, emptyGrids = 0,
		-- corners recovered by crossing two real face planes, and the ones where
		-- the planes were too near parallel to cross usefully
		corners = 0, unstableCorners = 0,
	}

	for part, g in pairs(localData.grids) do
		stats.parts += 1
		if g.fallback then stats.fallback += 1 else stats.block += 1 end
		local rings = ringsOfGrid(g, c, stats)
		if #rings == 0 then
			stats.emptyGrids += 1
			continue
		end
		-- Largest magnitude is the outer ring. Magnitude, never the sign: this
		-- project has been bitten by a handedness assumption before.
		local outer = rings[1]
		for _, r in ipairs(rings) do
			if math.abs(r.area) > math.abs(outer.area) then outer = r end
		end
		local region = {
			part = part,
			fallback = g.fallback,
			verts = outer.verts,
			area = math.abs(outer.area),
			cells = #g.cells,
			dead = #g.dead,
			holes = {},
		}
		for _, r in ipairs(rings) do
			if r ~= outer then
				region.holes[#region.holes + 1] = { verts = r.verts, area = math.abs(r.area) }
			end
		end
		stats.verts += #region.verts
		for _, h in ipairs(region.holes) do stats.verts += #h.verts end
		stats.regions += 1
		stats.holes += #region.holes
		regions[#regions + 1] = region
	end

	stats.seconds = os.clock() - t0
	return { regions = regions, stats = stats, config = c }
end

--------------------------------------------------------------------------

function Boundary.visualize(res: any, parent: Instance?)
	local root = parent or workspace
	local old = root:FindFirstChild("NavGen_Boundary")
	if old then old:Destroy() end
	local folder = Instance.new("Folder")
	folder.Name = "NavGen_Boundary"
	folder.Parent = root

	local function bar(a: Vector3, b: Vector3, colour: Color3, thick: number, into: Instance)
		local d = b - a
		if d.Magnitude < 1e-4 then return end
		local p = Instance.new("Part")
		p.Anchored = true; p.CanCollide = false; p.CanQuery = false; p.CanTouch = false
		p.Material = Enum.Material.Neon
		p.Color = colour
		p.Size = Vector3.new(thick, thick, d.Magnitude)
		p.CFrame = CFrame.lookAt(a + d * 0.5, b) + Vector3.new(0, 0.08, 0)
		p.Parent = into
	end

	-- Outer rings red, holes cyan, and a FALLBACK part's outline in amber --
	-- those are the world-aligned ones, the only place a staircased rim can
	-- still appear, so they should be identifiable without counting anything.
	for ri, r in ipairs(res.regions) do
		local sub = Instance.new("Folder")
		sub.Name = string.format("R%d_%s_c%d_h%d%s", ri, r.part.Name, r.cells, #r.holes,
			r.fallback and "_FALLBACK" or "")
		sub.Parent = folder
		local outerCol = r.fallback and Color3.fromRGB(255, 170, 60) or Color3.fromRGB(255, 80, 80)
		local v = r.verts
		for i = 1, #v do
			bar(v[i], v[i % #v + 1], outerCol, 0.18, sub)
		end
		for _, hr in ipairs(r.holes) do
			local hv = hr.verts
			for i = 1, #hv do
				bar(hv[i], hv[i % #hv + 1], Color3.fromRGB(80, 200, 255), 0.18, sub)
			end
		end
	end
	return folder
end

return Boundary
